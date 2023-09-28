import Icon from './svg/brain.svg';

function BrainIcon({ size }: { size: number }) {
  return (
    <div style={{ width: size }}>
      <img src={Icon} alt="brain icon" />
    </div>
  );
}

export default BrainIcon;
